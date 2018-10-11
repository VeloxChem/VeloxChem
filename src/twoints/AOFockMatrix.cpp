//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "AOFockMatrix.hpp"

#include <cmath>

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
        
        // spin restricted closed-shell Hatree-Fock
        
        if (dmtyp == denmat::rest)
        {
            _fockMatrices.push_back(CDenseMatrix(nrow, ncol));
            
            _fockTypes.push_back(fockmat::restjk);
            
            _scaleFactors.push_back(1.0);
            
            _idDensityMatrices.push_back(i);
        }
        
        // FIX ME: Add unrestricted open-shell Hatree-Fock
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

void
CAOFockMatrix::zero()
{
    for (size_t i = 0; i < _fockMatrices.size(); i++)
    {
        _fockMatrices[i].zero(); 
    }
}

int32_t
CAOFockMatrix::getNumberOfFockMatrices() const
{
    return static_cast<int32_t>(_fockMatrices.size());
}


int32_t
CAOFockMatrix::getNumberOfRows(const int32_t iFockMatrix) const
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        return _fockMatrices[iFockMatrix].getNumberOfRows();
    }
    
    return 0;
}

int32_t
CAOFockMatrix::getNumberOfColumns(const int32_t iFockMatrix) const
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        return _fockMatrices[iFockMatrix].getNumberOfColumns();
    }
    
    return 0;
}

int32_t
CAOFockMatrix::getNumberOfElements(const int32_t iFockMatrix) const
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        return _fockMatrices[iFockMatrix].getNumberOfElements();
    }
    
    return 0;
}

const double*
CAOFockMatrix::getFock(const int32_t iFockMatrix) const
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        return _fockMatrices[iFockMatrix].values();
    }
    
    return nullptr;
}

double*
CAOFockMatrix::getFock(const int32_t iFockMatrix)
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        return _fockMatrices[iFockMatrix].values();
    }
    
    return nullptr;
}

fockmat
CAOFockMatrix::getFockType(const int32_t iFockMatrix) const
{
    return _fockTypes[iFockMatrix];
}

double
CAOFockMatrix::getScaleFactor(const int32_t iFockMatrix) const
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        return _scaleFactors[iFockMatrix];
    }
    
    return 0.0;
}

int32_t
CAOFockMatrix::getDensityIdentifier(const int32_t iFockMatrix) const
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        return _idDensityMatrices[iFockMatrix];
    }
    
    return -1;
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
operator<<(      std::ostream&     output,
           const CAOFockMatrix& source)
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
