//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "AODensityMatrix.hpp"
#include "ErrorHandler.hpp"

CAODensityMatrix::CAODensityMatrix()

    : _denType(denmat::rest)
{
    
}

CAODensityMatrix::CAODensityMatrix(const std::vector<CDenseMatrix>& denMatrices,
                                   const denmat                     denType)

    : _denMatrices(denMatrices)

    , _denType(denType)
{
    if (denType == denmat::unrest)
    {
        errors::assertMsgCritical(
                denMatrices.size() % 2 == 0,
                "Odd number of matrices for unrestricted AO density!"
                );
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
CAODensityMatrix::getNumberOfDensityMatrices() const
{
    // restricted density matrix
    
    if (_denType == denmat::rest)
    {
        return static_cast<int32_t>(_denMatrices.size());
    }
    
    // unrestricted density matrix
    
    if (_denType == denmat::unrest)
    {
        return static_cast<int32_t>(_denMatrices.size()) / 2;
    }
    
    return 0;
}


int32_t
CAODensityMatrix::getNumberOfRows(const int32_t iDensityMatrix) const
{
    // restricted density matrix
    
    if (_denType == denmat::rest && iDensityMatrix < getNumberOfDensityMatrices())
    {
        return _denMatrices[iDensityMatrix].getNumberOfRows();
    }
    
    // unrestricted density matrix
    
    if (_denType == denmat::unrest && iDensityMatrix < getNumberOfDensityMatrices())
    {
        return _denMatrices[2 * iDensityMatrix].getNumberOfRows();
    }
    
    return 0;
}

int32_t
CAODensityMatrix::getNumberOfColumns(const int32_t iDensityMatrix) const
{
    // restricted density matrix
    
    if (_denType == denmat::rest && iDensityMatrix < getNumberOfDensityMatrices())
    {
        return _denMatrices[iDensityMatrix].getNumberOfColumns();
    }
    
    // unrestricted density matrix
    
    if (_denType == denmat::unrest && iDensityMatrix < getNumberOfDensityMatrices())
    {
        return _denMatrices[2 * iDensityMatrix].getNumberOfColumns();
    }
    
    return 0;
}

int32_t
CAODensityMatrix::getNumberOfElements(const int32_t iDensityMatrix) const
{
    // restricted density matrix
    
    if (_denType == denmat::rest && iDensityMatrix < getNumberOfDensityMatrices())
    {
        return _denMatrices[iDensityMatrix].getNumberOfElements();
    }
    
    // unrestricted density matrix
    
    if (_denType == denmat::unrest && iDensityMatrix < getNumberOfDensityMatrices())
    {
        return _denMatrices[2 * iDensityMatrix].getNumberOfElements();
    }
    
    return 0;
}

const double*
CAODensityMatrix::totalDensity(const int32_t iDensityMatrix) const
{
    if (_denType == denmat::rest && iDensityMatrix < getNumberOfDensityMatrices())
    {
        return _denMatrices[iDensityMatrix].values();
    }
    
    return nullptr;
}

const double*
CAODensityMatrix::alphaDensity(const int32_t iDensityMatrix) const
{
    if (_denType == denmat::unrest && iDensityMatrix < getNumberOfDensityMatrices())
    {
        return _denMatrices[2 * iDensityMatrix].values();
    }
    
    return nullptr;
}

const double*
CAODensityMatrix::betaDensity(const int32_t iDensityMatrix) const
{
    if (_denType == denmat::unrest && iDensityMatrix < getNumberOfDensityMatrices())
    {
        return _denMatrices[2 * iDensityMatrix + 1].values();
    }
    
    return nullptr;
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

std::ostream&
operator<<(      std::ostream&     output,
           const CAODensityMatrix& source)
{
    output << std::endl;
    
    output << "[CAODensityMatrix (Object):" << &source << "]" << std::endl;
    
    output << "_denType: " << to_string(source._denType); 
    
    output << "_denMatrices: " << std::endl;
    
    for (size_t i = 0; i < source._denMatrices.size(); i++)
    {
        output << "_denMatrices[" << i << "]: "<< std::endl;
        
        output << source._denMatrices[i] << std::endl;
    }
    
    return output;
}
